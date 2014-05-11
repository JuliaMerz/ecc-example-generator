function send(){
    var data = $("#ecc-form").serializeArray();
    $.post("ecc/", data, function(data, textStatus, jqXHR){
        $("#output").css("display", "block").empty().append(data);
    });
    return false;
}

